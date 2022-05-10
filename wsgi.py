#!/usr/bin/env python3

import sys
import site

site.addsitedir('/data/www/flask/genomedownload/venv/lib/python3.6/site-packages')
site.addsitedir('/data/www/flask/genomedownload/venv/lib/python3.6/site-packages/flask')
sys.path = list(filter(('/data/www/flask/genomedownload/').__ne__, sys.path))
sys.path.insert(0, '/data/www/flask/genomedownload/')



from app_dg import app as application

