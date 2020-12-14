# hyperspy-filtering

The application was created for the purpose of my engineer's thesis. The idea of the program was to make an easy and convinient tool for materials engineers that would help them analyse materials' microstructure images in an automated matter. The logic of the application is fully based on the code of Piotr Macioł, PhD. Eng.

## Installation and running

In order to make the application runnable, please follow the below steps:

1. Download the latest release of Hyperspy-bundle from https://github.com/hyperspy/hyperspy-bundle/releases
2. Run the installer that will install the distribution in the current directory.
3. Open the newly created folder, open WinPython Command Prompt and run the commands below:
   - `pip install seaborn`
   - `pip install django`
   - `pip install mpld3`
4. Move .bcf files you want to analyse to the "data" folder.
5. In WinPython Command Prompt, move to the application's directory (precisely the directory with Manage.py file) and run `python manage.py runserver` command.
