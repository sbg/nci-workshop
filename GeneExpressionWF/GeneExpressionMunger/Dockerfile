# Use the reference python container
FROM rfranklin/pythondev:latest

# Add the files from the GitHub repo [glob all (/*) files from each folder and put into homedir (/)]
RUN mkdir data
COPY munger.py /
COPY data/* data/

CMD ["/bin/bash"]