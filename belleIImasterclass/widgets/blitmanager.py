class BlitManager: 
    '''
    managing the blitting for the interactive widgets: only changed parts will be drawn in each iteration
    '''
    def __init__(self, canvas, artist, artist2=None):
        """copy from matplotlib website (blitting)"""
        self._canvas = canvas
        self._background = None
        self._artist = artist
        if artist2 is not None:
            self._artist2=artist2
            self._twoartists=True
        else:
            self._twoartists=False
        
        # grab the background on every draw
        self.cid = canvas.mpl_connect("draw_event", self.on_draw)

    def on_draw(self, event) -> None:
        canvas = self._canvas
        if event is not None:
            if event.canvas != canvas:
                raise RuntimeError
        self._background = canvas.copy_from_bbox(canvas.figure.bbox)
        self._draw_animated()


    def _draw_animated(self) -> None:
        fig = self._canvas.figure
        fig.draw_artist(self._artist)
        if self._twoartists:
            fig.draw_artist(self._artist2)

    def update(self) -> None:
        fig = self._canvas.figure
        if self._background is None:
            self.on_draw(None)
        else:
            self._canvas.restore_region(self._background)
            self._draw_animated()
            self._canvas.blit(fig.bbox)
        self._canvas.flush_events()


