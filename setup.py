import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="triangle_calculations_yes",
    version="0.0.2",
    author="LinQiuyun",
    author_email="haizhimenhao@outlook.com",
    description="Triangle calculations",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yesandnoandperhaps/Triangle",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=[
        'numpy', 
    ],
)

