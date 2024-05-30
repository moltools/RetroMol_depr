import requests
import zipfile
import os

url = "https://antismash.secondarymetabolites.org/upload/example/"

# Get all zip files in the directory at url

response = requests.get(url)
html = response.text

# Find all zip files
zips = []
for line in html.split("\n"):
    if ".zip" in line:
        zips.append(line.split('"')[1])

# report on zip files found
print("Found", len(zips), "zip files")
for i, z in enumerate(zips):
    print(i, z)

for zip_file in zips:

    output = zip_file
    foldername = zip_file.split(".")[0]

    response = requests.get(os.path.join(url, zip_file))

    with open(output, "wb") as f:
        f.write(response.content)

    # unzip
    with zipfile.ZipFile(output, 'r') as zip_ref:
        zip_ref.extractall(foldername)

    # delete old zip
    os.remove(output)

    print("Downloaded file to", output)