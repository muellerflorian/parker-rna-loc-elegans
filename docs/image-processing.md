# Image processing
Below we describe basic image processing operations, which are used in several workflows.

### Multi-channel conversion
All workflows require that each channel is stored as a separate z-stack. 
If images are stored as multi-channel z-stacks, they should be split in individual 
channels. This can usually be performed directly in the image acquisition software, or 
also with the freely available software
<a href="https://fiji.sc/"  target="_blank">**FIJI**</a>
as described below.

1.  Open image stack in FIJI, e.g. your `.dv`, `.lif`, or `.tif` images
0.  Split channels (FQ only supports mono-channel input images)
    a. From menu: `Image` > `Color` > `Split channels`
    a. Save each channels with a unique channel identifier, e.g. with a prefix `CH1_` or `DAPI_`.
0.  Depending on the particular workflow, images belonging you might need to store images belonging to one field
    of view in a separate folder.
