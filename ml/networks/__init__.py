# The code below ensures that all networks are automatically loaded and defined
# so that we can # iterate over the networks and automatically train/evaluate
# all of them. So the developer would not need to define and work with multiple
# notebooks and so on, they only need to implement the network

# Source - https://stackoverflow.com/a
# Posted by fortune_pickle
# Retrieved 2025-12-19, License - CC BY-SA 4.0

from inspect import isclass
from pkgutil import iter_modules
from pathlib import Path
from importlib import import_module

# iterate through the modules in the current package
package_dir = Path(__file__).resolve().parent
for (_, module_name, _) in iter_modules([package_dir]):

    if module_name.startswith("x"):
        continue

    # import the module and iterate through its attributes
    module = import_module(f"{__name__}.{module_name}")
    for attribute_name in dir(module):
        attribute = getattr(module, attribute_name)

        if isclass(attribute):
            # Add the class to this package's variables
            globals()[attribute_name] = attribute
