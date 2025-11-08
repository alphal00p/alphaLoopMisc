try:
    from symbolica import *
except ImportError:
    # User may not have selected symbolica features to be available
    pass

import re

floating_point_re = re.compile(r'((?:\d+(?:\.\d*)?|\.\d+)(?:`\d*)?)[eE]([+-]?\d+)')


def to_symbolica(expr: str, debug=False) -> Expression:
    if debug:
        print("Resulting string form of expression from Mathematica:\n%s"%expr)
    input_str = expr.replace(
        'MATHEMATICABACKTICK', '`').replace('"', 'MATHEMATICADOUBLEQUOTES').replace('*^','e')
    if debug:
        print("Expression sent for parsing with Symboolica:\n%s"%input_str)
    e = E(input_str, ParseMode.Mathematica)
    if debug:
        print("Resulting expression in symbolica:\n%s"%e.to_canonical_string())
    return E(input_str, ParseMode.Mathematica)


def to_mathematica_form(expr: str) -> str:
    processed_expr = expr.replace('MATHEMATICADOUBLEQUOTES', '"')
    return re.sub(floating_point_re, r'\1*^\2',  processed_expr)
