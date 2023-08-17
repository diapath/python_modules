# Python modules
Our Python modules

# ndpa.py

This modules can read and write NDPA annotations via a `ndpa.Annotation()` object.

One application might be to use [VALIS](https://github.com/MathOnco/valis) to transfer NDPA annotations from one NDPI WSI to another. The method for [transferring QuPath annotations](https://github.com/MathOnco/valis/issues/13) can easily be adapted to do so.

A point of detail, a NDPA file on its own does not store the pixel size (depending on the scanner lens used), nor the annotation physical offset. These are stored in the associated NDPI file and are read using the [tifffile module](https://pypi.org/project/tifffile/).

So, whether you use the `read_ndpa()` or `write_ndpa()` functions, the missing info will be read (if possible) from the associated .ndpi file, which is normally stored in the same directory. `/path/somefile.ndpi.ndpa` -> `/path/somefile.ndpi`

When drawing annotations in [NDP.View](https://www.hamamatsu.com/eu/en/product/life-science-and-medical-systems/digital-slide-scanner/U12388-01.html), we ask to put the class name (Tissue, Tumor, ...) either in the Title or in the Details field of an annotation. By default we look for the class name in the Details field (parameter `class_in_details=True`), and we make up a color based on the QuPath default class names and color scheme (parameter `import_colors=False`). These can be altered depending on intended use. An additional `distance_th` parameter can be passed on to simplify the annotation, but the method is very crude. If this is of interest, try to pass `distance_th=100` and `debug=True` to see the reduction in number of points.

For now, the recognised NDPA objects are: Freehand region, Freehand line, Circle, Rectangle, Ruler and Pin. Circles and rectangles are simply converted to polygons.

Also, since NDPA annotations are precious, to overwrite an existing NDPA annotation, you need to pass `overwrite=True` to the `write_ndpa()` otherwise an Exception is raised.

API may change, but for now:

## To read a NDPA file into a GeoJSON feature collection:

```python
import ndpa

fn = 'somefile.ndpi.ndpa'
annot = ndpa.Annotation()
fc = annot.read_ndpa(fn, class_in_details=False, import_colors=True)
```

## To write a GeoJSON feature collection into a NDPA file:

```python
import ndpa

fn = 'somefile.ndpi.ndpa'
annot = ndpa.Annotation()
annot.write_ndpa(fc, fn)
```

## An alternative method is to reference the NDPI file first:

```python
import ndpa

# Define a new annotation object and
# read its associated pixel size and offset from an NDPI file:
fn = 'somefile.ndpi'
annot = ndpa.Annotation(fn)

# To read the NDPA file associated with the NDPI file:
fc = annot.read_ndpa()

# To write a GeoJSON feature collection into the associated NDPA file:
annot.write_ndpa(fc)
```
