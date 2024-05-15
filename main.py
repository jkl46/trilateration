import re
import subprocess
import configparser
from kml import *

def getCoords(coordList):
    with open("demo.kml", "r") as f:
        s = f.read()

    res = re.findall("<coordinates>(.*)</coordinates>", s)
    coordList.extend([coord for group in res for coord in group.split(',')[::-1] if coord != '0'])

if __name__ == "__main__":
    coords = []
    getCoords(coords)

    subprocess.call(["cmake", "--build", "."])
    subprocess.call(["./main.exe", *coords])

    config = configparser.ConfigParser()
    config.read("trilaterate.ini")

    with open("results.kml", "w") as f:
        with kmlWrapper("trilateraiton results", f):
            for section in config.sections():
                with kmlFolder(section, f):
                    for key in config[section]:
                        _type = section.split(":")[-1]
                        name = key.replace("_", " ")
                        match _type:
                            case "point":
                                lat, lon = config[section][key].split(",")
                                marker = getPlaceMark(_type, name, lat, lon)

                            case "circle":
                                marker = getPlaceMark(_type, name, config[section][key])

                            case "line":
                                marker = ""
                                lat1, lon1, lat2, lon2 = [(x) for group in config[section][key].split(" ") for x in group.split(",")]
                                marker = getPlaceMark(_type, name, lat1, lon1, lat2, lon2)

                        f.write(marker)


