# 02413400

# importing relavant libraries:
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.image as mpimg
import csv # for writing

"""
Banksy Image Colour Changing triangles top and bottom
"""
# Provide the correct relative or absolute path to your image
image_path = '/Users/hedgehog/Desktop/MechEng_Degree/ME2_All/Computing_ME2/Python_ME2/NumericalMethods/NumericalMethods/img_Data/Flower.jpg'
img_matrix = mpimg.imread(image_path) # Load the image as a NumPy array
img_mod = img_matrix.copy()  # creates a writable version

# 1) Display the image:
plt.imshow(img_matrix)
plt.axis('off')  # Hide axes
#plt.show()

# 2) Render blue and red
# img_matrix is 333x500 matrix and rgb, top triangle only b, bottom only r
slope = (np.shape(img_matrix)[0])/(np.shape(img_matrix)[1])
# looping through each pixel
for i in range(0,(np.shape(img_matrix)[0])): # going down height
    for j in range(0,(np.shape(img_matrix)[1])): # going along width half way
        # Top triangle → BLUE (remove red and green)
        if i < slope * j and i < (np.shape(img_matrix)[0]) - slope * j:
            img_mod[i, j, 0] = 0  # red off
            img_mod[i, j, 1] = 0  # green off

        # Bottom triangle → RED (remove green and blue)
        elif i > slope * j and i > (np.shape(img_matrix)[0]) - slope * j:
            img_mod[i, j, 0] = 0  # green off
            img_mod[i, j, 1] = 0  # blue off

#plt.imsave('mattiesflower.jpg',img_mod)


"""
Stripy Image Colour Changing
"""
# Provide the correct relative or absolute path to your image
image_path = '/Users/hedgehog/Desktop/MechEng_Degree/ME2_All/Computing_ME2/Python_ME2/NumericalMethods/NumericalMethods/img_Data/Flower.jpg'
img_matrix = mpimg.imread(image_path) # Load the image as a NumPy array
img_mod = img_matrix.copy()  # creates a writable version

# 1) Display the image:
plt.imshow(img_matrix)
plt.axis('off')  # Hide axes
#plt.show()

# 2) Render blue and red
# img_matrix is 333x500 matrix and rgb, top triangle only b, bottom only r
slope = (np.shape(img_matrix)[0])/(np.shape(img_matrix)[1])
# looping through each pixel
for i in range(0,(np.shape(img_matrix)[0])): # going down height
    for j in range(0,(np.shape(img_matrix)[1])): # going along width half way
        # Top triangle → BLUE (remove red and green)
        if i < slope * j and i < (np.shape(img_matrix)[0]) - slope * j:
            img_mod[i, j, 0] = 0  # red off
            img_mod[i, j, 1] = 0  # green off

        # Bottom triangle → RED (remove green and blue)
        elif i > slope * j and i > (np.shape(img_matrix)[0]) - slope * j:
            img_mod[i, j, 0] = 0  # green off
            img_mod[i, j, 1] = 0  # blue off

plt.imsave('mattiesflower.jpg',img_mod)