# experiments with the Python Image Library (PIL)
# free from:  http://www.pythonware.com/products/pil/index.htm
# create 128x128 (max size) thumbnails of all JPEG images in the working folder

import glob
import Image


def convert_png_image_to_thumbnail(str):
    nstr = str.replace(".png", ".thumbnail.png");
    im = Image.open(str);
    im.thumbnail((128, 128), Image.ANTIALIAS);
    im.save(nstr, "png");

#get image files

for infile in glob.glob("*.png"):
    convert_png_image_to_thumbnail(infile);
