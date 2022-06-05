FROM python:3.10.4
COPY ./rosspy
RUN pip install -r ./requirements.txt