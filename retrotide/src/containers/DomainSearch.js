import React from 'react';
import ModuleBuilder from '../components/ModuleBuilder'

class DomainSearch extends React.Component {

  constructor(props) {
    super(props);
    this.state = {
      PKSDomainList: {
        CAL: {domainName: 'CAL', present: false},
        KS:  {domainName: 'KS', present: false},
        AT:  {domainName: 'AT', present: false},
        DH:  {domainName: 'DH', present: false},
        ER:  {domainName: 'ER', present: false},
        KR:  {domainName: 'KR', present: false},
        ACP: {domainName: 'ACP', present: false}, 
        TE:  {domainName: 'TE', present: false},
      },
      DomainList: 'PKS',
    }
  }

  // select PKS or NRBS
  // list of all modulebuilders
  // add module button
  // delete module button
  // submit button

  render() {
    return (
      <div className='DomainSearch'>
        <h3>Construct Modules</h3> 
        <ModuleBuilder domainList={this.state.PKSDomainList} />
      </div>
    )
  }
};

export default DomainSearch;