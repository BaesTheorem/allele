"""Read consumer DNA exports and produce an annotated report, entirely offline."""

from .model import Call, Sample
from .parsers import parse
from .promethease import parse_report

__version__ = "0.1.0"
__all__ = ["Call", "Sample", "parse", "parse_report", "__version__"]
