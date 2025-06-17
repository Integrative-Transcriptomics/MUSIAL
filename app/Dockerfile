# syntax=docker/dockerfile:1

FROM python:3.10-slim-buster

WORKDIR /app

ENV JAVA_HOME=/opt/java/openjdk
COPY --from=eclipse-temurin:17-jre $JAVA_HOME $JAVA_HOME
ENV PATH="${JAVA_HOME}/bin:${PATH}"

COPY requirements.txt requirements.txt
RUN pip3 install -r requirements.txt

COPY . .

ENV PYTHONUNBUFFERED=1

EXPOSE 5001
CMD ["gunicorn","--workers", "4","--timeout","600","--graceful-timeout", "300", "--bind", "0.0.0.0:5001","musialweb:app"]