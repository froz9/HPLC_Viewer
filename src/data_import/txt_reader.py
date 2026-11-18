"""Parser for tab-delimited and general text chromatographic files (.txt)."""

from __future__ import annotations

import io
from pathlib import Path
from typing import Optional, Union

from src.data_import.csv_reader import (
    DelimitedParseResult,
    parse_delimited_chromatogram,
)


def read_txt_chromatogram(
    source: Union[str, Path, bytes, io.BytesIO],
    file_name: Optional[str] = None,
    raise_on_error: bool = False,
) -> DelimitedParseResult:
    """Read a tab-delimited or general delimited TXT chromatographic file (.txt)."""
    return parse_delimited_chromatogram(
        source=source,
        file_name=file_name,
        preferred_delimiter="\t",
        raise_on_error=raise_on_error,
    )
