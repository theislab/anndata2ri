from __future__ import annotations

from rpy2.robjects import conversion, numpy2ri


converter = conversion.Converter('original scipy conversion', template=numpy2ri.converter)
