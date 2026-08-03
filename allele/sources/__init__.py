"""Annotation sources.

Three, deliberately different in kind:

- ClinVar    clinical significance, curated by submitters and expert panels
- CPIC       pharmacogenomic prescribing guidance, gene level
- GWAS       statistical trait associations from published studies
- SNPedia    subjective interest grading from a wiki, via your own report

They disagree, which is why the merge in `annotate.py` keeps all three rather
than picking one.
"""

from .base import Annotation, AnnotationSource, SourceInfo
from .clinvar import ClinVar
from .cpic import Cpic
from .gwas import GwasCatalog
from .snpedia import SNPedia

__all__ = [
    "Annotation",
    "AnnotationSource",
    "SourceInfo",
    "ClinVar",
    "Cpic",
    "GwasCatalog",
    "SNPedia",
]
