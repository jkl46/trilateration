kmlPrepend = """<?xml version="1.0" encoding="UTF-8"?>

<kml xmlns="http://www.opengis.net/kml/2.2" xmlns:gx="http://www.google.com/kml/ext/2.2" xmlns:kml="http://www.opengis.net/kml/2.2" xmlns:atom="http://www.w3.org/2005/Atom">
<Document>
	<name>Main point.kml</name>
	<StyleMap id="m_ylw-pushpin">
		<Pair>
			<key>normal</key>
			<styleUrl>#s_ylw-pushpin</styleUrl>
		</Pair>
		<Pair>
			<key>highlight</key>
			<styleUrl>#s_ylw-pushpin_hl</styleUrl>
		</Pair>
	</StyleMap>
	<Style id="s_ylw-pushpin">
		<IconStyle>
			<scale>1.1</scale>
			<Icon>
				<href>http://maps.google.com/mapfiles/kml/pushpin/ylw-pushpin.png</href>
			</Icon>
			<hotSpot x="20" y="2" xunits="pixels" yunits="pixels"/>
		</IconStyle>
	</Style>
	<Style id="s_ylw-pushpin_hl">
		<IconStyle>
			<scale>1.3</scale>
			<Icon>
				<href>http://maps.google.com/mapfiles/kml/pushpin/ylw-pushpin.png</href>
			</Icon>
			<hotSpot x="20" y="2" xunits="pixels" yunits="pixels"/>
		</IconStyle>
	</Style>
		"""
kmlAppend = """
</Document>
</kml>
	"""


def getPlaceMark(type, *kw):
	match type:
		case "point":
			name, lat, lon = kw
			body = """
				<Placemark>
					<name>%s</name>
					<visibility>1</visibility>
					<styleUrl>#m_ylw-pushpin</styleUrl>
					<Point>
						<gx:drawOrder>1</gx:drawOrder>
						<coordinates>%s, %s</coordinates>
					</Point>
				</Placemark>
			""" % (name, lon, lat)
		case "circle":
			name, coordString= kw
			body = """
				<Placemark>
				<name>%s</name>
				<visibility>1</visibility>
				<styleUrl>#inline</styleUrl>
				<LineString>
					<tessellate>1</tessellate>
					<coordinates>
						%s
					</coordinates>
				</LineString>
			</Placemark>	
			""" % (name, coordString)
		case "line":
			name, lat1, lon1, lat2, lon2 = kw
			body = """
	<Placemark>
			<name>%s</name>
			<styleUrl>#inline</styleUrl>
			<LineString>
				<tessellate>1</tessellate>
				<coordinates>
				%s,%s,0 %s,%s,0
				</coordinates>
			</LineString>
		</Placemark>
			""" % (name, lon1, lat1, lon2, lat2)

	return body

class kmlFolder():
	def __init__(self, folderName, fileHandle):
		self.folderName = folderName
		self.fileHandle = fileHandle

	def __enter__(self):
		text = """<Folder>
		<name>%s</name>
		<open>1</open>\n""" % self.folderName
		self.fileHandle.write(text)

	def __exit__(self, *kw):
		text = "</Folder>\n"
		self.fileHandle.write(text)

class kmlWrapper():
	def __init__(self, fileName, fileHandle):
		self.fileName = fileName
		self.fileHandle = fileHandle

	def __enter__(self):
		text = """<?xml version="1.0" encoding="UTF-8"?>

<kml xmlns="http://www.opengis.net/kml/2.2" xmlns:gx="http://www.google.com/kml/ext/2.2" xmlns:kml="http://www.opengis.net/kml/2.2" xmlns:atom="http://www.w3.org/2005/Atom">
<Document>
	<name>%s</name>
	<StyleMap id="m_ylw-pushpin">
		<Pair>
			<key>normal</key>
			<styleUrl>#s_ylw-pushpin</styleUrl>
		</Pair>
		<Pair>
			<key>highlight</key>
			<styleUrl>#s_ylw-pushpin_hl</styleUrl>
		</Pair>
	</StyleMap>
	<Style id="s_ylw-pushpin">
		<IconStyle>
			<scale>1.1</scale>
			<Icon>
				<href>http://maps.google.com/mapfiles/kml/pushpin/ylw-pushpin.png</href>
			</Icon>
			<hotSpot x="20" y="2" xunits="pixels" yunits="pixels"/>
		</IconStyle>
	</Style>
	<Style id="s_ylw-pushpin_hl">
		<IconStyle>
			<scale>1.3</scale>
			<Icon>
				<href>http://maps.google.com/mapfiles/kml/pushpin/ylw-pushpin.png</href>
			</Icon>
			<hotSpot x="20" y="2" xunits="pixels" yunits="pixels"/>
		</IconStyle>
	</Style>
		""" % self.fileName
		self.fileHandle.write(text)

	def __exit__(self, *kw):
		text = """
</Document>
</kml>
	"""
		self.fileHandle.write(text)
    


