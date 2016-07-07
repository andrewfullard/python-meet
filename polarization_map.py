from PIL import Image 


image_05_10k_0 = Image.open('/home/manisha/Research/diff_temp/bow_05_0_10k/bow_05_0_10k_0.jpeg')
image_05_10k_90 = Image.open('/home/manisha/Research/diff_temp/bow_05_0_10k/bow_05_0_10k_90.jpeg')
image_05_20k_0 = Image.open('/home/manisha/Research/diff_temp/bow_15_0_10k/bow_15_0_10k_0.jpeg')
image_05_20k_90 = Image.open('/home/manisha/Research/diff_temp/bow_15_0_10k/bow_15_0_10k_90.jpeg')

colorbar = image_05_10k_0.crop((0, 10, 10, 10))

croppedIm_05_10k_0 = image_05_10k_0.crop((100, 150, 550, 650))
croppedIm_05_10k_90 = image_05_10k_90.crop((100, 150, 550, 650))
croppedIm_05_20k_0 = image_05_20k_0.crop((100, 150, 550, 650))
croppedIm_05_20k_90 = image_05_20k_90.crop((100, 150, 550, 650))

width,height = croppedIm_05_10k_0.size
croppedIm_05_10k_0 = croppedIm_05_10k_0.resize((width/2,height/2))
croppedIm_05_10k_90 = croppedIm_05_10k_90.resize((width/2,height/2))
croppedIm_05_20k_0 = croppedIm_05_20k_0.resize((width/2,height/2))
croppedIm_05_20k_90 = croppedIm_05_20k_90.resize((width/2,height/2))
colorbar = colorbar.resize((width/2,height/2))
width1,height1 = croppedIm_05_10k_0.size

polmap = Image.new('RGB', (width, height))
#polmap = image_05_10k_0.copy()
polmap.paste(croppedIm_05_10k_0,(0,0))
polmap.paste(croppedIm_05_10k_90,(0,height1))
polmap.paste(croppedIm_05_20k_0,(width1,0))
polmap.paste(croppedIm_05_20k_90,(width1,height1))

polmap.save('/home/manisha/Research/polmap_10k_0.eps', format='eps', dpi=1000)
colorbar.save('/home/manisha/Research/colorbar.eps', format='eps', dpi=100)
polmap.show()
