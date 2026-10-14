import gzip
import subprocess
import tempfile
from io import BytesIO
from pathlib import Path

import pytest

from malva.xopen import xopen


def has_bgzip():
    """Check if bgzip is available in the system."""
    try:
        subprocess.run(["bgzip", "--version"], capture_output=True, check=True)
        return True
    except (subprocess.SubprocessError, FileNotFoundError):
        return False


class TestXopen:
    @pytest.fixture
    def temp_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    def test_file_like_objects(self):
        content = b"File-like\nObject\n"

        # BytesIO
        bytes_io = BytesIO(content)
        with xopen(bytes_io, "rb") as f:
            assert f.read() == content

    @pytest.mark.skipif(not has_bgzip(), reason="bgzip not available")
    def test_bgzip_compression(self, temp_dir):
        input_file = temp_dir / "test.fastq"
        output_file = temp_dir / "test.fastq.gz"
        content = "@read1\nACGT\n+\nIIII\n"

        with open(input_file, "w") as f:
            f.write(content)

        subprocess.run(
            ["bgzip", "-c", str(input_file)], stdout=open(str(output_file), "wb"), check=True
        )

        with xopen(output_file, "r") as f:
            assert f.read() == content

    def test_pipe_handling(self):
        # Test reading from stdin
        with pytest.raises(ValueError):
            with xopen("-", "r") as f:
                f.read()

    def test_read_text_file(self, temp_dir):
        text_file = temp_dir / "test.txt"
        content = "Hello\nWorld\n"

        with open(text_file, "w") as f:
            f.write(content)

        with xopen(text_file, "r") as f:
            assert f.read() == content

    def test_read_gzip_file(self, temp_dir):
        gz_file = temp_dir / "test.gz"
        content = "Compressed\nContent\n"

        with gzip.open(gz_file, "wt") as f:
            f.write(content)

        with xopen(gz_file, "r") as f:
            assert f.read() == content

    def test_append_modes(self, temp_dir):
        text_file = temp_dir / "append.txt"

        with xopen(text_file, "w") as f:
            f.write("First\n")

        with xopen(text_file, "a") as f:
            f.write("Second\n")

        with open(text_file, "r") as f:
            assert f.read() == "First\nSecond\n"

    def test_binary_modes(self, temp_dir):
        binary_file = temp_dir / "binary.dat"
        content = b"Binary\nContent\n"

        with xopen(binary_file, "wb") as f:
            f.write(content)

        with xopen(binary_file, "rb") as f:
            assert f.read() == content

    def test_error_handling(self, temp_dir):
        # Test invalid mode
        with pytest.raises(ValueError):
            xopen(temp_dir / "test.txt", "x")  # ty: ignore[no-matching-overload]

        # Test invalid compression format
        with pytest.raises(ValueError):
            xopen(temp_dir / "test.txt", "r", format="invalid")

        # Test non-existent file
        with pytest.raises(FileNotFoundError):
            xopen(temp_dir / "nonexistent.txt", "r")

    def test_bzip2_roundtrip(self, temp_dir):
        path = temp_dir / "test.bz2"
        with xopen(path, "w") as f:
            f.write("bzip2\ncontent\n")
        with xopen(path, "r") as f:
            assert f.read() == "bzip2\ncontent\n"

    def test_xz_roundtrip(self, temp_dir):
        path = temp_dir / "test.xz"
        with xopen(path, "w") as f:
            f.write("xz\ncontent\n")
        with xopen(path, "r") as f:
            assert f.read() == "xz\ncontent\n"

    def test_format_override_without_extension(self, temp_dir):
        path = temp_dir / "noext"
        with xopen(path, "w", format="gz") as f:
            f.write("override\n")
        with xopen(path, "rb", format="gz") as f:
            assert f.read() == b"override\n"

    def test_threads_zero_uses_python_backend(self, temp_dir):
        path = temp_dir / "threads.gz"
        with xopen(path, "w", threads=0) as f:
            f.write("threads\n")
        with xopen(path, "r", threads=0) as f:
            assert f.read() == "threads\n"

    def test_content_detection_without_extension(self, temp_dir):
        path = temp_dir / "detect.dat"
        with gzip.open(path, "wt") as f:
            f.write("detected\n")
        with xopen(path, "r") as f:
            assert f.read() == "detected\n"

    def test_compresslevel_out_of_range(self, temp_dir):
        with pytest.raises(ValueError):
            xopen(temp_dir / "x.gz", "w", compresslevel=99, threads=0)

    def test_unsupported_format(self, temp_dir):
        with pytest.raises(ValueError):
            xopen(temp_dir / "x.txt", "r", format="zip")


if __name__ == "__main__":
    pytest.main([__file__])
