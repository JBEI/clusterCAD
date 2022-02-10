import React from 'react';
import Button from '../components/Button';

class DomainSearch extends React.Component {

  constructor(props) {
    super(props);
    this.state = {modules: ''}
    console.log(this.state);
  }

  render() {
    return (
      <div className='DomainSearch'>
        <h3>Construct Module</h3>
        <div className="inputSelector">
          <Button> A
          </Button>
        </div>
      </div>
    )
  }
};

export default DomainSearch;