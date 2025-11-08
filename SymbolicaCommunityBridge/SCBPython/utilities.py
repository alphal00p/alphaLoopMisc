try:
    from symbolica import *
except ImportError:
    # User may not have selected symbolica features to be available
    pass

def to_symbolica(expr: str) -> Expression:

    input_str = expr.replace(
        'MATHEMATICABACKTICK', '`').replace('"', 'MATHEMATICADOUBLEQUOTES')
    #print(input_str)
    return E(input_str, ParseMode.Mathematica)


def to_mathematica_form(expr: str) -> str:
    return expr.replace('MATHEMATICADOUBLEQUOTES', '"')
