#   Copyright 2023 DIAPath - CMMI Gosselies Belgium (diapath@cmmi.be)
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.import numpy as np

import xml.etree.ElementTree as ET
from xml.dom import minidom
import tifffile
import json
import uuid
import os

class Feature(dict):
    def __init__(self, name=None, geometry_type='Polygon', classification=None, c=-1, z=0, t=0):
        self['type'] = 'Feature'
        self['id'] = str(uuid.uuid4())
        self['geometry'] = {'type':geometry_type,'coordinates':[]}
        self['properties'] = {
            'object_type': 'annotation',
            'isLocked': True
        }
        if name is not None:
            self['properties']['name'] = name
        if classification is not None:
            self['properties']['classification'] = classification

        # Plane information (the plane to which the annotation belongs)
        self['plane'] = {
            "c": c,
            "z": z,
            "t": t
        }

    @property
    def coordinates(self):
        return self['geometry']['coordinates']

    @coordinates.setter
    def coordinates(self,coordinates):
        self['geometry']['coordinates'] = coordinates

    @property
    def size(self):
        return len(self['geometry']['coordinates'])

class Classification(dict):
    def __init__(self, name, color=None):
        """name: the classification name (e.g. Tumor)
        if the name is recognised and no color is given, the one from qupath will be given
        color: the classification color as a hex number, but becomes / saved as colorRGB
        """

        #in case no color was passed
        if color is None:
            #These are the QuPath default class names and colors
            colormap = {"Tumor":0xc80000, "Tissue":0x00c800, "Stroma":0x96c896, "Immune cells":0xa05aa0,
                    "Necrosis":0x323232, "Region*":0x0000b4, "Ignore*":0xb4b4b4, "Positive":0xfa3e3e,
                    "Negative":0x7070e1, "Other":0xffc800, "Clear":0x000000}

            if name in colormap.keys():
                color = colormap[name]
            else:
                # a default red color...
                color=0xc80000

        self['name'] = name
        self.color = color

    @property
    def color(self):
        return self['colorRGB'] & 0xffffff
       
    # a setter function
    @color.setter
    def color(self, c):
        if type(c) == tuple or type(c) == list:
            #[r,g,b]
            c = c[0] & 0xff << 4 + c[1] & 0xff << 2 + c[2] & 0xff

        self['colorRGB'] = c - 0x1000000


class Annotation:
    px_nm,py_nm = None, None
    xoffset_nm,yoffset_nm = None,None
    debug = False
    ndpi_path = None

    def __init__(self, ndpi_path=None, debug=False):
        """An NDPI file is required to get the pixel size, image centre and image offset.
        However, these parameters can be defined directly via set_params()
        """

        if ndpi_path is not None:
            self.ndpi_path = ndpi_path
            if not ndpi_path.lower().endswith(".ndpi"):
                raise ValueError('Not an NDPI file...')

            self.get_ndpi_params(ndpi_path, debug)

    def set_params(self, px_nm=1,py_nm=1, xoffset_nm=0,yoffset_nm=0):
        self.px_nm = px_nm
        self.py_nm = py_nm
        self.xoffset_nm = xoffset_nm
        self.yoffset_nm = yoffset_nm

    def get_ndpa_path(self, ndpa_path=None):
        if ndpa_path is None:
            if self.ndpi_path is None:
                raise ValueError('No reference NDPI file defined...')
            else:
                ndpa_path = self.ndpi_path + ".ndpa"

        if not ndpa_path.lower().endswith('.ndpa'):
            raise ValueError('Not an NDPA file...')

        return ndpa_path

    def read_ndpa(self, ndpa_path=None, class_in_details=True, import_colors=False, distance_th=0., debug=None):
        """Here we read a NDPA file

        Notes:
        Here the issue is that NDPA does not define holes in polygons. We can define a name "Clear" for holes
        NDPA come with two strings, the title and the details, either of which can be used to define a class.

        :param ndpa_path: The path of the NDPA file to be read.
        :param class_in_details: The 'details' field contains the name of the class
        :param import_colors: NDPA colours are read, otherwise use QuPath colors associated with class name
        :param distance_th: Polygon simplification when > 0 (in pixel size units, a good value is 100 nm)
        :param debug: output extra debug information
        """

        if debug is None: debug = self.debug

        ndpa_path = self.get_ndpa_path(ndpa_path)
        if not os.path.exists(ndpa_path):
            raise ValueError('NDPA annotation file not found.')

        if self.px_nm is None:
            ndpi_path = ndpa_path[:-5]
            if not os.path.exists(ndpi_path):
                raise ValueError('Reference NDPI image not found.')

            # Get the parameters from the associated NDPI file
            self.get_ndpi_params(ndpi_path, debug)

        #The feature collection we return
        fc = {
            "type": "FeatureCollection",
            "features": []
        }
        
        root = ET.parse(ndpa_path).getroot()
        nodes = root.findall('ndpviewstate')
        n_nodes = len(nodes)
        for i,ndpviewstate_node in enumerate(nodes):
            annotation_name = ndpviewstate_node.find('title').text
            annotation_details = ndpviewstate_node.find('details').text
            if annotation_details is None:
                annotation_details = ''

            if class_in_details:
                class_name = annotation_details.strip().capitalize()
            else:
                class_name = annotation_name

            #FIXME not sure what these are for...
            annotation_x = ndpviewstate_node.find('x').text
            annotation_y = ndpviewstate_node.find('y').text
            annotation_z = ndpviewstate_node.find('z').text

            annotation_node = ndpviewstate_node.find('annotation')
            annotation_type = annotation_node.attrib['type'].upper()
            annotation_displayname = annotation_node.attrib['displayname']

            if import_colors:
                annotation_color = int(annotation_node.attrib['color'][1:],16) - 0x1000000
            else:
                annotation_color = None

            if debug:
                print(f"ndpviewstate node {i+1}/{n_nodes}: Name is {annotation_name}, class is {class_name}, color is {annotation_color}")

            #by default, we all the annotations are closed (polygons)
            annotation_isclosed = 1

            #TODO Add any other annotation types not already dealt with. These are handled:
            # <annotation type="freehand" displayname="AnnotateFreehandLine" color="#000000">
            # <annotation type='freehand' displayname='AnnotateFreehand' color='#c80000'>
            # <annotation type="freehand" displayname="AnnotateRectangle" color="#00ff00">
            # <annotation type="circle" displayname="AnnotateCircle" color="#ff0000">
            # <annotation type="pin" displayname="AnnotatePin" color="#ffffff">

            if annotation_type == "FREEHAND":
                pointlist_node = annotation_node.find('pointlist')

                #Points are common to all the nodes.
                points = []
                for point_node in pointlist_node.findall('point'):
                    x, y = float(point_node.find('x').text), float( point_node.find('y').text)  # physical position
                    x,y = self.ndpacoords_to_coords(x,y)
                    points.append((int(x), int(y)))

                if distance_th > 0:
                    points = self.reduce_polygon(points, distance_th=distance_th, debug=debug)

                annotation_isclosed = int(annotation_node.find('closed').text)

                if annotation_isclosed  == 1:
                    geometry_type = 'Polygon'
                    if str(points[0]) != str(points[-1]):
                        points.append(points[0])
                    points = [points]

                else:
                    geometry_type = 'LineString'

                #specialtype can indicate a rectangle...
                if "specialtype" in annotation_node.attrib and annotation_node.attrib["specialtype"] == "rectangle":
                    print("We seem to have a rectangle... We don't really care")

            elif annotation_type == "CIRCLE":
                geometry_type = 'Polygon'
                x = int(annotation_node.find('x').text)
                y = int(annotation_node.find('y').text)
                x,y = self.ndpacoords_to_coords(x,y)
                r = int(int(annotation_node.find('radius').text)/self.px_nm)

                if debug: print(f"Circle: x={x}, y={y}, r={r}")
                t = np.linspace(-np.pi, np.pi, 50)
                points = [np.array([r*np.cos(t)+x, r*np.sin(t)+y],int).transpose().tolist()]

            elif annotation_type == "PIN":
                geometry_type = 'Point'
                x = int(annotation_node.find('x').text)
                y = int(annotation_node.find('y').text)
                x,y = self.ndpacoords_to_coords(x,y)
                # Here we have a single point https://en.wikipedia.org/wiki/GeoJSON
                points = [x,y]

            elif annotationType == "LINEARMEASURE":
                geometry_type = 'LineString'
                x1 = int(annotation_node.find('x1').text)
                y1 = int(annotation_node.find('y1').text)
                x1,y1 = self.ndpacoords_to_coords(x1,y1)
                x2 = int(annotation_node.find('x2').text)
                y2 = int(annotation_node.find('y2').text)
                x2,y2 = self.ndpacoords_to_coords(x2,y2)
                points = [[x1,y1],[x2,y2]]

            else:
                print(f"Problem !!! Annotation type {annotation_type} not handled right now...")
                continue

            # Here the annotation name serves as the class name (can be Tissue, etc...)
            c = Classification(class_name, annotation_color)
            f = Feature(annotation_name, geometry_type, classification=c, z=annotation_z)
            f.coordinates = points

            fc['features'].append(f)

        return fc

    
    def write_ndpa(self, fc, ndpa_path=None, overwrite=False, hole_name='Clear', debug=None):
        """
        Write the given GEOJSON feature collection to a NDPA file

        :param fc: GEOJSON feature collection to save as NDPA file
        :param ndpa_path: Path to the output NDPA file.
        :param overwrite: By default, do not overwrite existing NDPA files.
        :param hole_name: The class name for holes in polygons (default is 'Clear').
        :param debug: display debug messages.
        """

        if debug is None: debug = self.debug

        if fc is None:
            raise ValueError('No feature collection to save!')

        ndpa_path = self.get_ndpa_path(ndpa_path)
        if os.path.exists(ndpa_path) and overwrite==False:
            raise ValueError('NDPA annotation file exists. Will not overwrite unless explicitely instructed!')

        if fc['type'] != 'FeatureCollection':
            print("Only JSON feature collections as input...")
            return

        # Get the parameters from the associated NDPI file
        if self.px_nm is None:
            ndpi_path = ndpa_path[:-5]
            if not os.path.exists(ndpi_path):
                raise ValueError('Reference NDPI image not found.')

            # Get the parameters from the associated NDPI file
            self.get_ndpi_params(ndpi_path, debug)
        
        root = ET.Element('annotations')
        for i, f in enumerate(fc['features']):
            if not f['type'] == 'Feature':
                print(f"Element {i+1} is not a feature...")
                continue

            # here we check the type of annotation
            g = f['geometry']
            annotation_type = g['type']
            if debug: print(f"Iterating feature {i+1}/{len(fc['features'])} ({annotation_type})")

            if  annotation_type == 'Polygon' or annotation_type == 'MultiPolygon':
                ndpa_type = 'freehand'
                ndpa_displayname = 'AnnotateFreehand'
                ndpa_closed = 1
                ndpa_coords = g['coordinates']
                if annotation_type == 'Polygon':
                    ndpa_coords = [ndpa_coords]

            elif annotation_type == 'LineString' or annotation_type == 'MultiLineString':
                ndpa_type = 'freehand'
                ndpa_displayname = 'AnnotateFreehandLine'
                ndpa_closed = 0
                ndpa_coords = g['coordinates']
                if annotation_type == 'MultiLineString':
                    ndpa_coords = [ndpa_coords]
                if annotation_type == 'LineString':
                    ndpa_coords = [[ndpa_coords]]

            elif annotation_type == 'Point' or annotation_type == 'MultiPoint':
                ndpa_type = 'pin'
                ndpa_displayname = 'AnnotatePin'
                ndpa_coords = g['coordinates']
                if annotation_type == 'MultiPoint':
                    ndpa_coords = [ndpa_coords]
                if annotation_type == 'Point':
                    ndpa_coords = [[ndpa_coords]]
            else:
                print(f"FIXME! {annotation_type} not supported yet!")
                continue

            # FIXME Here we can get the z, but is it needed?... I haven't seen any variation in <ndpviewstate><z>, it's always 0 even if changing plane before drawing a new annotation.
            #if 'plane' in g.keys():
            #    annotation_z = g['plane']['z'] if 'z' in g['plane'].keys() else 0
            annotation_z = 0

            # the annotation name (different from the classname)
            annotation_name = ''
            class_name = ''
            c = 0xc80000

            if 'properties' in f.keys():
                if 'name' in f['properties'].keys():
                    anntotation_name = f['properties']['name']
                if 'classification' in f['properties'].keys():
                    if 'name' in f['properties']['classification'].keys():
                        class_name = f['properties']['classification']['name']
                    if 'color' in f['properties']['classification'].keys():
                        c = f['properties']['classification']['color']
                    elif 'colorRGB' in f['properties']['classification'].keys():
                        c = f['properties']['classification']['colorRGB']
            ndpa_color = self.color_to_hexstring(c)

            #FIXME is it safe to leave x and y as 0 0? This means no need to store xcenter and ycenter
            dic = {
                'title': annotation_name,
                'details': class_name,
            }

            for j, coords_list in enumerate(ndpa_coords):
                if debug: print(f"Iterating coord_list {j+1}/{len(ndpa_coords)}")
                for k, coords in enumerate(coords_list):
                    ndpviewstate_node = self.make_ndpviewstate_node(root, i, **dic)

                    # make the annotation node
                    annotation_node = ET.SubElement(ndpviewstate_node, 'annotation',
                        {'type': ndpa_type, 'displayname':ndpa_displayname, 'color':ndpa_color})
                    ET.SubElement(annotation_node, 'measuretype').text = '0'

                    #Here we have different behaviours depending on the class type...
                    if ndpa_type == 'freehand':
                        ET.SubElement(annotation_node, 'closed').text = str(ndpa_closed)
                        if k > 0:
                            # These are holes, we need to change the details to "Clear" and the color to black
                            annotation_node.set('color', '#000000')
                            ET.SubElement(ndpviewstate_node, 'details').text = hole_name

                        #if closed, do not repeat first and last point...
                        if ndpa_closed == True and str(coords[0]) == str(coords[-1]):
                            coords = coords[:-1]

                        #here we add the pointlist
                        pointlist_node = ET.SubElement(annotation_node, 'pointlist')
                        for point in coords:
                            point_node = ET.SubElement(pointlist_node, 'point')
                            x, y = self.coords_to_ndpacoords(*point)
                            ET.SubElement(point_node, 'x').text = str(x)
                            ET.SubElement(point_node, 'y').text = str(y)

                    elif ndpa_type == 'pin':
                        x, y = self.coords_to_ndpacoords(*coords)
                        ET.SubElement(annotation_node, 'x').text = str(x)
                        ET.SubElement(annotation_node, 'y').text = str(y)
                        ET.SubElement(annotation_node, 'icon').text = 'placered'
                        ET.SubElement(annotation_node, 'stricon').text = 'iconplacered'

        xmlstr = minidom.parseString(ET.tostring(root, method="html", encoding='UTF-8')).toprettyxml(indent="   ")
        with open(ndpa_path, "w") as f:
            f.write(xmlstr)

    def reduce_polygon(self, polygon, angle_th=0, distance_th=0, debug=None):
        # code lifted from StackOverflow: https://stackoverflow.com/a/69026086

        if debug is None: debug=self.debug

        n_before = len(polygon)
        angle_th_rad = np.deg2rad(angle_th)
        points_removed = [0]
        while len(points_removed):
            points_removed = list()
            for i in range(0, len(polygon)-2, 2):
                v01 = np.array(polygon[i-1]) - np.array(polygon[i])
                v12 = np.array(polygon[i]) - np.array(polygon[i+1])

                d01 = np.linalg.norm(v01)
                d12 = np.linalg.norm(v12)
                if d01 < distance_th and d12 < distance_th:
                    points_removed.append(i)
                    continue

                val = np.sum(v01*v12) / (d01 * d12)
                if val <= 1:
                    angle = np.arccos(val)
                    if angle < angle_th_rad:
                        points_removed.append(i)

            polygon = np.delete(polygon, points_removed, axis=0)

        if debug:
            print(f"reduce_polygon before / after: {n_before} -> {len(polygon)}")

        return polygon.tolist()

    def coords_to_ndpacoords(self, x, y):
        """Convert pixel coordinates to physical NDPA coordinates"""

        x, y = x * self.px_nm, y * self.py_nm
        x, y = x + self.xoffset_nm, y + self.yoffset_nm
        return int(x), int(y)

    def ndpacoords_to_coords(self, x, y):
        """Convert physical NDPA coordinates to pixel coordinates"""

        x, y = x - self.xoffset_nm, y - self.yoffset_nm
        x, y = x / self.px_nm, y / self.py_nm
        return int(x), int(y)

    def get_ndpi_params(self, ndpi_path, debug=None):
        """Get pixel size and image offset from an NDPI file.
        Unfortunately these are not contained in the NDPA metadata,
        so an NDPA file is always associated to an NDPI file.
        """

        if debug is None: debug=self.debug

        tif = tifffile.TiffFile(ndpi_path)

        ru = tif.pages[0].tags['ResolutionUnit']
        if not 'CENTIMETER' in str(ru):
            # I'm sure we could be more cleverer about this... CENTIMETER -> 1e-2 in xyres, MICROMETER... etc etc.. but is it needed?
            # That's why I'm raising a ValueError, let me know if this ever becomes an issue
            raise ValueError('ResolutionUnit should be \'CENTIMETER\'!')
        
        #Resolution unit is centimeter
        xres = tif.pages[0].tags['XResolution'].value[0]
        px_nm = (1e-2 / xres) * 1e9

        #Resolution unit is centimeter
        yres = tif.pages[0].tags['YResolution'].value[0]
        py_nm = (1e-2 / yres) * 1e9

        #FIXME Each <ndpviewstate> seems to need a <lens> and <x>,<y>,<z> parameters
        #Previously I used the image xcenter and ycenter for x and y, but 0,0 (and 0 for z) seem to work
        #just as well. If not, here's xcenter and ycenter:
        #self.xcenter = tif.pages[0].tags['ImageWidth'].value // 2
        #self.ycenter = tif.pages[0].tags['ImageLength'].value // 2

        xcenter_nm = (tif.pages[0].tags['ImageWidth'].value / 2) * px_nm
        ycenter_nm = (tif.pages[0].tags['ImageLength'].value / 2) * py_nm

        #Offset between NDPA coordinates and 0 based pixel coordinates. Offset units in nm
        xoffset_nm = tif.pages[0].ndpi_tags['XOffsetFromSlideCenter'] - xcenter_nm 
        yoffset_nm = tif.pages[0].ndpi_tags['YOffsetFromSlideCenter'] - ycenter_nm

        self.set_params(px_nm,py_nm,xoffset_nm,yoffset_nm)

        if debug:
            print("Pixel size (nm):",px_nm, py_nm)
            print("Image center (nm):", xcenter_nm,ycenter_nm)
            print("Combined Offset (nm):",xoffset_nm,yoffset_nm)
            
        tif.close()

    def color_to_hexstring(self, c, debug=None):
        """color can either be an integer or a (r,g,b) tuple or list"""

        if debug is None: debug=self.debug

        if debug: print("Color input:", c)

        if type(c) == tuple or type(c) == list:
            #[r,g,b]
            c = ((c[0] & 0xff) << 16) + ((c[1] & 0xff) << 8) + (c[2] & 0xff)
        
        c = '#'+hex(c & 0xffffff)[2:]
        if debug: print("Color output:", c)

        return c

    def make_ndpviewstate_node(self, root, i, **kwargs):
        node = ET.SubElement(root, 'ndpviewstate', {'id': str(i + 1)})

        #these would be the default values for a ndpviewstate children
        dic_defaults = {
            'title': '',
            'details': '',
            'coordformat':'nanometers',
            'lens':0.445623,
            'x':0,
            'y':0,
            'z':0,
            'showtitle':0,
            'showhistogram':0,
            'showlineprofile':0,
        }

        for k,v in dic_defaults.items():
            if k in kwargs.keys(): v = kwargs[k]
            ET.SubElement(node, k).text = str(v)

        return node

