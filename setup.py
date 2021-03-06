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
    long_description_content_type="text/markdown",
    author="Omochice",
    author_email="h20ms419@hirosaki-u.ac.jp",
    url="https://github.com/MizLab/mizlab-tools",
    license="MIT",
    install_requires=[
        "biopython",
        "requests",
        "nptyping",
        "numpy",
    ],
    packages=find_packages(exclude=(
        "tests",
        "docs",
        "source",
    )),
    tests_require=[
        "pytest",
        "pytest-runner",
        "pytest-cov",
    ],
    entry_points={
        "console_scripts": [
            "fetch_gbk = mizlab_tools.fetch_gbk:main",
            "fetch_taxon = mizlab_tools.fetch_taxon:main",
            "calculate_coordinates = mizlab_tools.calculate_coordinates:main",
            "calculate_weights = mizlab_tools.calculate_weights:main",
            "weight_to_table= mizlab_tools.weight_to_table:main"
        ]
    },
    python_requires=">=3.6.1",
)
