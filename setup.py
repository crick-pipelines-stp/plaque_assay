from setuptools import setup


def read_requirements():
    with open("requirements.txt") as f:
        return [i.strip() for i in f.readlines()]


setup(
    name="plaque_assay",
    version="0.3",
    author="Scott Warchal",
    author_email="scott.warchal@crick.ac.uk",
    license="MIT",
    packages=["plaque_assay", "plaque_assay.titration"],
    python_requires=">=3.6",
    install_requires=read_requirements(),
    tests_require="pytest",
    zip_safe=True,
)
