
import re

# Function to format a single monomial term
def format_monomial(monomial):
    match = re.match(r"([+-]?)(\d+)?(?:\*)?([a-zA-Z0-9^\*]+)?(?:/(\d+))?", monomial)
    if not match:
        return None

    sign = match.group(1) or ''
    numerator = match.group(2) or '1'
    monomial_part = match.group(3) or ''
    denominator = match.group(4) or '1'
    
    fraction_part = f"{numerator}/{denominator}" if denominator != '1' else numerator
    fraction_part = fraction_part if numerator != '1' or denominator != '1' else ''
    return f"{sign}{fraction_part}{('*' + monomial_part) if monomial_part and fraction_part else monomial_part}"

# Function to format a polynomial
def format_polynomial(poly):
    poly=poly.replace('**','^').replace(' ','')
    if not poly.startswith('+') and not poly.startswith('-'):
        poly = '+' + poly
    terms = re.split(r"([+-])", poly)[1:]
    terms = [terms[i] + terms[i+1] for i in range(0, len(terms), 2)]
    formatted_terms = [format_monomial(term) for term in terms]
    return ''.join(formatted_terms).lstrip('+')

# Test function for monomials
def test_format_monomial():
    test_cases = [
        ("22434*x*y^2/12", "22434/12*x*y^2"),
        ("-22434*x*y^2/12", "-22434/12*x*y^2"),
        ("+22434*x*y^2/12", "+22434/12*x*y^2"),
        ("24/3", "24/3"),
        ("77*x*y^2", "77*x*y^2"),
        ("x*y^2/44", "1/44*x*y^2"),
        ("x^18*z", "x^18*z"),
    ]
    for case, expected in test_cases:
        result = format_monomial(case)
        assert result == expected, f"Failed for {case}, got {result}, expected {expected}"
    print("All monomial tests passed.")

# Test function for polynomials
def test_format_polynomial():
    test_cases = [
        ("a/3+k1x^2/7", "1/3*a+1/7*k1x^2"),
        ("-a/3+k1x^2/7", "-1/3*a+1/7*k1x^2"),
        ("a/3-k1x^2/7", "1/3*a-1/7*k1x^2"),
        ("22434*x*y^2/12+24/3-77*x*y^2+x*y^2/44+x^18*z", "22434/12*x*y^2+24/3-77*x*y^2+1/44*x*y^2+x^18*z"),
        ("-42*x*y^123+y/33", "-42*x*y^123+1/33*y")
    ]
    for case, expected in test_cases:
        result = format_polynomial(case)
        assert result == expected, f"Failed for {case}, got {result}, expected {expected}"
    print("All polynomial tests passed.")

if __name__ == "__main__":
    test_format_monomial()
    test_format_polynomial()
