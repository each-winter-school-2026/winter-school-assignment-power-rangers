<!-- Open the visual viewer using Ctrl + shift + v -->
















































# EACH Winter School 2026 exercise
## Setting up

When the codespace was launched all dependencies should automatically have been installed. 

Start the EWOKS server by typing the following in the terminal:
```bash
python app/manage.py runserver

```

Which should result in something similar to this

```
Watching for file changes with StatReloader
Performing system checks...

System check identified no issues (0 silenced).
January 15, 2026 - 13:38:59
Django version 5.2.8, using settings 'EWOKS.settings'
Starting development server at http://127.0.0.1:8000/
Quit the server with CTRL-BREAK.

WARNING: This is a development server. Do not use it in a production setting. Use a production WSGI or ASGI server instead.
For more information on production servers see: https://docs.djangoproject.com/en/5.2/howto/deployment/
```

The server can be stopped by hitting `ctrl + c` in the terminal where the server is running.

>When using the visual interface the backend logs are shown here. These are primarily for debugging and can usually be ignored.

The server can then be opened by clicking on `http://127.0.0.1:8000/` within the terminal.

This should open a separate browser tab which should display the visual interface of EWOKS. 

### Troubleshooting

On startup the codespace should have automatically installed all dependencies. 
If this is not the case they can be installed by running the following two lines in the terminal

```bash
pip install -r requirements.txt

export PYTHONPATH=${PWD}
```

---

## Assignment

EWOKS is a visual *In silico* proteomics modeling pipeline designed to model analytical methods within code. 
A visual interface is provided so that researchers without programming knowledge can utilize it to predict and verify the outcome of experiments.

In this assignment you will be documenting and modeling an analytical method of your choice. You can select an analytical method with your team within the [google spreadsheet](https://edu.nl/c9yxn).
If your favorite method is taken and none of the other methods interest you, you can add your own to the sheet. When you have selected a method with your team ensure that you mark it as taken within the sheet.

### Step 1: Identify a suitable model for *in silico* prediction

Find a simple model that can simulate your chosen analytical method on-the-fly. Consider the following:

**Performance Requirements:**
- The model should be fast to compute (no more than a few seconds for an entire proteome of ~20,000 proteins)
- Sequence-based models are preferred as they are easier to implement and typically perform faster than structure-based models

**Model Options:**
- **Sequence-only models**: Try to find models that work purely from amino acid sequence (preferred)
- **Structure-based models**: If sequence-only simulation isn't feasible, AlphaFold-predicted structures are available from UniProt (explore this if you have time)
- **Existing AI/deep learning models**: If you find suitable pre-trained models created by others, document these methods for potential integration

### Step 2: Write a short description for the visual interface

Your description should be short and concise (max 100 words). It should clearly convey the function of your module.

> **Tip:** GitHub Codespaces has GitHub Copilot built in. This is a powerful AI tool which you may use to your advantage.

### Step 3: Implement the model in Python if possible

A demonstration on how to implement a module can be found in [HOW_TO_IMPLEMENT_A_MODULE.md](HOW_TO_IMPLEMENT_A_MODULE.md)

### Step 4: Find and document the parameters of commercially available kits

- Research the parameters of commercially available kits related to your chosen analytical method.
- Extensively document these commercial kits. Include all information that is required to simulate this kit *in silico*.
- If possible implement these values as an option the user can choose within your model.

### Step 5: Prepare a presentation (Thursday morning)

Prepare a short presentation about your chosen analytical method.

At a minimum you should include:

- A brief explanation on how your analytical method works

---

## Need Help?

If you encounter any issues or have questions during the assignment:
- Ask your team members for assistance
- Consult the documentation files in this repository
- Reach out to the instructors for guidance

Good luck!