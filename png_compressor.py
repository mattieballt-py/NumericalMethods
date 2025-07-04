from PIL import Image

def compress_png(input_path, output_path, optimize=True, reduce_colors=False):
    with Image.open(input_path) as img:
        if reduce_colors:
            img = img.convert('P', palette=Image.ADAPTIVE)

        img.save(output_path, optimize=optimize)

        original_size = round(os.path.getsize(input_path) / 1024, 2)
        compressed_size = round(os.path.getsize(output_path) / 1024, 2)

        print(f"✅ Compression complete: {original_size} KB ➜ {compressed_size} KB")

# Example usage
import os

input_file = "/Users/hedgehog/Desktop/MechEng_Degree/ME2_All/Computing_ME2/Python_ME2/NumericalMethods/NumericalMethods/img_Data/Screenshot 2025-07-04 at 18.28.00 1.png"
output_file = "your_image_compressed.png"
compress_png(input_file, output_file, optimize=True, reduce_colors=True)
