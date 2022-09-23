import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="HavNegpy",
    version="0.1",
    author="Mohamed Kolmangadi",
    author_email="mohamed.kolmangadi@gmail.com",
    description="Python package to analyze dielectric data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=["HavNegpy"],
    install_requires=["numpy", "scipy", "mplcursors", "matplotlib"]
)
