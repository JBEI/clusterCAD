import React from 'react';
import { connect } from 'react-redux';
import { updateMolecule } from '../redux/actions/actions';
import Button from './Button';
import ChemDoodleWrapper from './ChemDoodleWrapper';

// This tab contains the canvas for drawing a target molecule and submitting it to backend
// It does not store the molecule information locally but makes a async call to create a
// job in the store and await results

const mapDispatchToProps = dispatch => {
  return {
    updateMolecule: molecule => dispatch(updateMolecule(molecule)),
    dispatch,
  }
};

class MoleculeDiagram extends React.Component {

  render() {
    return (
      <div className='MoleculeDiagram form'>
        <ChemDoodleWrapper />
        <p>Use the ChemDoodle canvas provided to draw the diagram of your target molecule</p>
        {/*<Button className='submitMolecule' disabled={false} onClick={send_test} >Submit</Button>*/}
        <Button className='submitMolecule' disabled={false} onClick={this.props.updateMolecule} >Test</Button>
      </div>
    )
  }
}


export default connect(null, mapDispatchToProps)(MoleculeDiagram);