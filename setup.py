import setuptools

setuptools.setup(
    name="RetroMol",
    version="0.0.1",
    author="David Meijer",
    author_email="david.meijer@wur.nl",
    install_requires=[],
    package_dir={"": "src"},
    packages=[
        "retromol",
        "retromol_sequencing"
    ],
    python_requires=">=3.10",
    entry_points={"console_scripts": ["retromol = main:main"]}
)