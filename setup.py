from setuptools import setup, find_packages

VERSION = "1.0.0"
DESCRIPTION = "Calculate shannon entropy of biological sequences"

setup(
    name="seqshannon",
    version=VERSION,
    author="Wojciech Dobrychłop",
    author_email="wojciech.dobrychlop@gmail.com",
    url="https://github.com/WojDob/SeqShannon",
    description=DESCRIPTION,
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    packages=find_packages(),
    install_requires=["biopython"],
    keywords=["python", "bioinformatics", "shannon entropy", "sequence analysis"],
    entry_points={
        "console_scripts": [
            "seqshannon = seqshannon.seqshannon:main",
        ],
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
    ],
)
