# RetroTideUI

# Project Structure

The Retrotide project lives inside the ClusterCAD project's folder and is entirely self-contained in this folder; it does not share any files or overlap with ClusterCAD. The intention is that a front-end developer can work on Retrotide webapp code without disrupting any of the code of the larger ClusterCAD project. Both projects run together in a Docker container as described below, but from the end user's perspective this is all one website and project. They also share a git repository at `https://github.com/JBEI/clusterCAD`. Updates to the Retrotide frontend are currently located in a branch called `react-app` and its offshoots. There is a branch called `domain-isolation` which contains the most up-to-date code for domain search only and has had other components deleted - this branch will need to be folded back into the main `react-app` branch at some point carefully without deleting all the other components. Note that `domain-isolation`'s stylesheets have been heavily edited to mimic the look of the ClusterCAD project - this will also have to be reconciled with the rest of retrotide if we opt not to make the websites identical.

The future plan for this project is that ClusterCAD frontend will be gradually moved over to React, possibly template by template, until the entire project is running as one single service on the same host.

# Retrotide in Docker Container

The entire ClusterCAD project is inside a single Docker container, detailed in the `docker-compose.yml` file one directory up. This project, Retrotide, is run as a service inside that container, named `react`. The service lives at host `:3000`. NB that ClusterCAD itself runs on the `:8000` port. These are the internal Docker ports but are exposed to the user for purposes of development. 
We are using proxies to allow the two projects to communicate - this is detailed in the `nginx` folder in the main `clusterCAD` project folder. The utility file is located in `sites-enabled/django_project`. This file specifies that any request coming in to a url prepended with `/retrotide` should be proxied over to the Retrotide project. Retrotide's `static` folder files must also be prepended with this string to force it to be separated from the Django project's `static` folder. Because React doesn't allow us to rename the `static` folder, all request urls must be prepended this way. This is enforced via a file in the `Retrotide` folder called `.env` which holds the environment variables and forces all file names in the `public` folder to be preparended with `/retrotide`. 
The `/retrotide` string is also added automatically to internal urls by the Redux-router found in the `App.js` file. We have set `/retrotide` as the `basename`. To navigate within the project, use Redux-router's `Link` tag which will automatically add the `/retrotide` for you. To leave the project use a regular `a` anchor tag. This means that to make a request to the API, we simply return control to nginx by using an `a` tag with no `/retrotide` string. 

# Static files and images

Images and other such files should be placed in `src -> images` and React will sort them into the appropriate static location. Icons are taken from Remix Icon 2.5.0 at https://remixicon.com/

# Running Retrotide Locally

In the root project directory, do `docker compose up` to bring up the entire Docker container, then navigate to `localhost:8000/retrotide` to view Retrotide running in your browser. Note that Docker hot-reload is currently not working, so updates to js and css files won't be seen on page refresh. You'll have to stop and remove the container in order to see changes.

To start Retrotide as a standalone React webapp, do `npm start` and navigate to `localhost:3000`. NB that this will only run the frontend with no backend availability and many functions won't work as expected. For developing in css or making other frontend changes this is faster than waiting for the entire Docker container, but all changes should be double-checked in the container before publishing. 

# Using The App

Retrotide has a few different search tools.

# Domain Search

Domain Search allows you to construct a PKS protein domain by domain. Navigate to `localhost:8000/retrotide/domainSearch`. There should be two modules already present, the loading and terminating modules. These cannot be edited or removed in the current version. To add more extending modules, click the "Add Modules +" button. To remove these extending modules, click the "X" button inside the module. Within extending modules, groups of compatible domains can be added or removed by clicking the "+" in the list of domains. Only one group can be added at a time. Some domains also have editable options - these domains appear in green. To edit options, click the domain name which will open the options modal. The currently selected option appears in dark green. Click another option and close the modal to change the selected option. To submit the finished PKS, click the "Submit" button at the top of the page. 