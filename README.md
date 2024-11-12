# OpenFOAM and ParaView Tips

---

### OpenFOAM Tips

- **See the Generated Mesh (works only outside `linuxsh`):**
  1. Open the correct directory.
  2. Run:
     ```bash
     blockMesh
     paraFoam &    # "&" allows the terminal to continue working
     ```
  3. In ParaView, select **'surface with edges'**. Note that the result shown is just the first iteration/timestep.

- **Save Log File (using `icoFoam` and running it in the background):**
  1. Generate the mesh.
  2. Run:
     ```bash
     icoFoam > log &
     ```

---

### ParaView Tips

- **Record Steps:**  
  Navigate to `Tools` > `Start Trace`. After stopping the trace, save it as a **macro**. This can be opened directly later.

- **Quick Select Filter:**  
  Use `Ctrl + Space`, type the filter name, and press `Enter`.

---

### Post-Processing

- **List Available Functions:**  
  Run:
  ```bash
  postProcess -list
