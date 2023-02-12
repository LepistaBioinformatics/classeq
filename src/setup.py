#!/usr/bin/python3
# -*- coding: UTF-8 -*-

__version__ = "0.1.0"

import sys

try:
    from setuptools import find_packages, setup
except ImportError:
    sys.exit(
        "We need the Python library setuptools to be installed. "
        + "Try running: python -m ensurepip"
    )


with open("../README.md", "r") as readme_file:
    README = readme_file.read()


if __name__ == "__main__":
    setup(
        author="Biotrop - Bioinformatic Team",
        author_email="bioinfo@biotrop.com.br",
        name="qs",
        version=__version__,
        description="Biotrop Quorum-Sensing",
        long_description=README,
        long_description_content_type="text/markdown",
        keywords=["Biotrop", "Bioinformatics"],
        package_data={"": ["assets/*.tsv.gz"]},
        packages=find_packages(),
        classifiers=[
            "Programming Language :: Python :: 3.8",
            "Operating System :: OS Independent",
        ],
        python_requires=">=3.8",
        setup_requires=["wheel"],
        entry_points={
            "console_scripts": [
                "qs=qs.ports.cli.main:qs_cmd",
                "qs-api=qs.ports.api.main:app",
            ],
        },
    )
