import React from 'react';
import {clusterCADSeqSearch} from '../redux/middleware/api';

const canvasDimensions = {
  HEIGHT: 500,
  WIDTH: 800,
}

class ChemDoodleWrapper extends React.Component {

  constructor(props) {
    super(props);
    this.state = {
      errorMessage: "",
      chemDoodleCanvas: {}, // JS object corresponding to sketcher
      molecule: "",
    }
  }

  componentDidMount() {
    let newCanvas = new window.ChemDoodle.SketcherCanvas(
      'chemDoodleCanvas',
      canvasDimensions.WIDTH,
      canvasDimensions.HEIGHT,
      {
        useServices:false,
        oneMolecule:true
      }
    );
    // we need to use a callback to get the result of setState as it's asynchronous
    this.setState({chemDoodleCanvas: newCanvas}, () => {console.log(this.state.chemDoodleCanvas)});
  }

  render() {
    return ( 
      <div id="chemDoodleWrapper">
        <canvas id="chemDoodleCanvas" width={canvasDimensions.WIDTH + "px"} height={canvasDimensions.HEIGHT + "px"} style={{ border: '1px solid black' }}>
          This is a message that shows if the browser doesn't support HTML5 canvas, which all modern browsers should support now.
        </canvas>
        <button onClick={() => {
          let mol = this.state.chemDoodleCanvas.getMolecule();
          let formattedMol = window.ChemDoodle.writeMOL(mol);
          this.setState({
              molecule: formattedMol,
            }, () => {clusterCADSeqSearch(this.state.molecule, "")}); // is it necessary to put molecule in state?
          }}>TEST</button>
        <p><small>
          {"Molecule sketcher is used courtesy of "}
          <a href="http://web.ChemDoodle.com" rel="noopener noreferrer" target="_blank">ChemDoodle</a>
          {" under version 3 of the GNU GPL license."}
        </small></p>
      </div>
    )
  }
}

export default ChemDoodleWrapper;