import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="solsys",
    version="0.0.1",
    author="Behrouz Safari",
    author_email="behrouz.safari@gmail.com",
    description="A python package for Solar System calculations",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/behrouzz/solsys",
    license="MIT",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    packages=["solsys"],
    include_package_data=True,
    install_requires=["numpy", "scipy", "pandas", "matplotlib", "requests", "jplephem", "hypatie"],
    python_requires='>=3.4',
)
