import re
from decimal import Decimal
from pprint import pprint

# Read the data from the file (Assuming the file is in the same directory as the script)
with open("out.ms", "r") as f:
    data = f.read()

# Regular expression replacement
# The pattern looks for integers separated by a space, a '/', another space, '2^', and then another integer.
# The replace function formats these in the desired Decimal structure.
def replace_decimal(match):
    int_num = match.group(1)
    base = match.group(2)
    exp = match.group(3)
    return f"Decimal({int_num} / {base}**{exp})"

modified_data = re.sub(r"(-?\d+) / (\d+)\^(\d+)", replace_decimal, data.replace(':',''))

# Evaluation
# Since we're using eval on data from an untrusted source, this could be risky in a different context.
# However, as we're using this for a specific scientific computation, the risk is assumed to be low.
evaluated_data = eval(modified_data)

# Pretty Printing
pprint(evaluated_data)
