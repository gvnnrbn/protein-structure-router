# Run app
uvicorn src.main:app --reload

# docker
docker build -t my-tapo-app .
docker run --rm -p 8000:8000 protein-router-test
