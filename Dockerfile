FROM golang
WORKDIR /app
RUN go get
ENTRYPOINT ["go", "run", "."]