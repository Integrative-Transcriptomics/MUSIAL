from musialweb import app
from flask import render_template


@app.route("/")
def get_template_home():
    return render_template("home.html")


@app.route("/results")
def get_template_results():
    return render_template("results.html")


@app.route("/about")
def get_template_about():
    return render_template("about.html")