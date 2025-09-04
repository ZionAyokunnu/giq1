from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="giq2",
    version="1.0.0",
    author="Zion (Jaron's lab)",
    supervisor="Kamil, Sasha, Arif, Sam",
    author_email="za7@sanger.ac.uk",
    description="A gene inversion algorithm for analyzing chromosomal rearrangements",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ZionAyokunnu/giq2",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 2 - Beta",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific :: Bio-Informatics",
    ],
    python_requires=">=3.8",
    install_requires=[
        "pandas>=1.3.0",
        "numpy>=1.20.0",
    ],
    extras_require={
        "dev": [
            "pytest>=6.0",
            "black>=21.0",
            "flake8>=3.8",
        ],
    },
    entry_points={
        "console_scripts": [
            "giq2=giq2.main:main",
        ],
    },
)
