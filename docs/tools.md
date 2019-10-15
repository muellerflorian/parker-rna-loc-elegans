

# Tools

Here we list the tools that are used in the different analysis workflows.

## FISH-quant: RNA detection

[**FISH-quant**](https://bitbucket.org/muellerflorian/fish_quant/) is a Matlab toolbox to
localize RNAs in 3D from smFISH images. FISH-quant is often used to perform the RNA
detection. Eventually, we would like to implement the main FISH-quant features in ImJoy.

<img src="https://raw.githubusercontent.com/muellerflorian/rna-loc/master/docs/img/fq-screenshot.png" width="600px"></img>


## ImJoy
[**ImJoy**](https://imjoy.io/docs/#/) is image processing platform with an easy to use interface powered by a Python engine running in the background. ImJoy plays a
central role in most analysis workflows.

<img src="https://raw.githubusercontent.com/muellerflorian/rna-loc/master/docs/img/imjoy-screenshot.png" width="600px"></img>

### Working with ImJoy
We provide links to install the different ImJoy plugins in dedicated **ImJoy workspaces**.
Workspaces can be selected from little puzzle symbol in the upper left part of the
interface.

Most plugins require the **ImJoy Plugin Engine**, to perform computations in
Python. You will need to **install** it only once, but **launch** it each time
you work with ImJoy. For more information for how to install and use the pluging engine, please consult the [ImJoy documentation](https://imjoy.io/docs/#/user-manual?id=python-engine).

Once installed, ImJoy remembers the workspaces and plugins and you simply have to open the web app and select the appropriate workspace [https://imjoy.io/#/app](https://imjoy.io/#/app)

If you press on the installation link, the ImJoy web app will open and display a dialog asking if you want to install the specified plugin. To confirm, press the `install` button.

<img  src="https://raw.githubusercontent.com/muellerflorian/rna-loc/master/docs/img/install_plugin.png" width="600px"></img>
