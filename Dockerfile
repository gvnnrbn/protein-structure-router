# Stage 1: Extract TAPO files from the original image
FROM pdoviet/tapo:v1.1.3-alpha.2 AS tapo_source

# Stage 2: Modern Python + Java runtime
FROM python:3.11-slim

# Install OpenJDK + libraries Swing needs for headless font rendering (TREks uses Swing)
RUN apt-get update && apt-get install -y --no-install-recommends \
    default-jre-headless \
    libharfbuzz0b \
    libfreetype6 \
    libfontconfig1 \
    wget \
    && rm -rf /var/lib/apt/lists/*

# Copy the TAPO application tree from stage 1
COPY --from=tapo_source /home/sgeadmin/save/BioApps/tapo-v1.1.3 /opt/tapo

# Copy the work directory: contains DSSP bank and local PDB cache
COPY --from=tapo_source /home/sgeadmin/work /home/sgeadmin/work

# JAXB was removed from Java 9+; BioJava 4.1.1 (bundled with TAPO) needs it.
# Drop the JARs into dependencies/ — already on the classpath via dependencies/* wildcard.
RUN wget -q "https://repo1.maven.org/maven2/javax/xml/bind/jaxb-api/2.3.1/jaxb-api-2.3.1.jar" \
        -O /opt/tapo/dependencies/jaxb-api-2.3.1.jar \
    && wget -q "https://repo1.maven.org/maven2/com/sun/xml/bind/jaxb-impl/2.3.9/jaxb-impl-2.3.9.jar" \
        -O /opt/tapo/dependencies/jaxb-impl-2.3.9.jar \
    && wget -q "https://repo1.maven.org/maven2/com/sun/activation/javax.activation/1.2.0/javax.activation-1.2.0.jar" \
        -O /opt/tapo/dependencies/javax.activation-1.2.0.jar

# Create writable directories TAPO expects
RUN mkdir -p /home/sgeadmin/work/tmp /opt/tapo/tmp /opt/tapo/log /opt/tapo/output \
    && chmod -R 777 /home/sgeadmin/work/tmp /opt/tapo/tmp /opt/tapo/log /opt/tapo/output /opt/tapo/data

# Set up the Python app
WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY src/ ./src/

ENV TAPO_DIR=/opt/tapo
ENV PORT=8000

CMD ["sh", "-c", "uvicorn src.main:app --host 0.0.0.0 --port ${PORT}"]