import os
from asyncio import sleep
from jupyter_server.base.handlers import JupyterHandler
import tornado

class XarrayLeafletHandler(JupyterHandler):

    def set_default_headers(self):
        self.set_header('Content-Type', 'image/png')
        self.set_header('Cache-Control', 'no-store, no-cache, must-revalidate, max-age=0')

    @tornado.web.authenticated
    async def get(self, path):
        path = '/' + path
        path_done = path[:-4] + '.done'
        delete = False
        timeout = False
        dt = 0.1
        t = 0
        while True:
            if os.path.exists(path_done):
                with open(path_done) as f:
                    txt = f.read()
                if txt.startswith('keep'):
                    break
                elif txt.startswith('delete'):
                    delete = True
                    break
            await sleep(dt)
            t += dt
            # don't wait more than 10 seconds
            # Leaflet can make a request while dragging but Python is not
            # triggered unless the mouse button is released, so the file
            # will never be written
            if t > 10:
                timeout = True
                break
            if t > 1:
                # we started polling every 0.1s, but after 1s we only poll
                # every second as we are not very reactive anyway!
                dt = 1
        if timeout:
            self.finish()
            return
        # serve the tile
        with open(path, 'rb') as f:
            tile_png = f.read()
        self.set_header('Cache-Control', 'no-store, no-cache, must-revalidate, max-age=0')
        self.finish(tile_png)
        if delete:
            os.remove(path_done)
            os.remove(path)
