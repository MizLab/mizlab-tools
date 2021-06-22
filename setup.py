from pathlib import Path

from setuptools import find_packages, setup


def load_readme() -> str:
    prj_root = Path(__file__).resolve().parent
    with open(str(prj_root / "README.md")) as f:
        return f.read()


setup(
    name="mizlab_tools",
    version="0.0.1",
    description="Some codes used in research.",
    long_description=load_readme(),
    author="Omochice",
    author_email="h20ms419@hirosaki-u.ac.jp",
    url="",
    license="MIT",
    install_requires=["biopython", "requests"],
    packages=find_packages(exclude=("tests", "docs")),
    tests_require=["pytest", "pytest-runner"],
    entry_points={
        "console_scripts": [
            "fetch_gbk = mizlab_tools.fetch_gbk:main",
            "fetch_taxon = mizlab_tools.fetch_taxon:main"
        ]
    },
)
