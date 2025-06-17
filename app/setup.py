from setuptools import setup

setup(
    name="musialweb",
    version="2.0",
    packages=["musialweb"],
    include_package_data=True,
    install_requires=[
        "flask==3.1.1",
		"Flask-Session==0.8.0",
		"werkzeug==3.1.3",
		"cachelib==0.13.0",
		"python-dotenv==1.1.0",
		"biopython==1.85",
		"pandas==2.3.0",
		"numpy==1.24",
		"scipy==1.15.3",
		"gunicorn==23.0.0",
		"scikit-learn==1.7.0",
		"igraph==0.11.9",
    ],
)
