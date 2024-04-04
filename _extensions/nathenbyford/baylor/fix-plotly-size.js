window.FixPlotlySize = {
  id: 'FixPlotlySize',
  init: function(Reveal) {
    Reveal.addEventListener( 'slidechanged', function( event ) {
      let curSlide = Reveal.getCurrentSlide();
      require(['plotly'], function(Plotly) {
        let plotlyDivs = curSlide.querySelectorAll('.js-plotly-plot');
        for (let i = 0; i < plotlyDivs.length; i++) {
          Plotly.Plots.resize(plotlyDivs[i]);
        }
      });
    });
  }
}