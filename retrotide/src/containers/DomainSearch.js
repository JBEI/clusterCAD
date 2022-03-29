import React from 'react';
import Button from '../components/Button';
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
      LoadingPKSList: {
        KS:  {domainName: 'KS', present: false},
        AT:  {domainName: 'AT', present: true},
        ACP: {domainName: 'ACP', present: true},       
      },
      TerminatingPKSList: {
        KS:  {domainName: 'KS', present: false},
        AT:  {domainName: 'AT', present: false},
        DH:  {domainName: 'DH', present: false},
        ER:  {domainName: 'ER', present: false},
        KR:  {domainName: 'KR', present: false},
        ACP: {domainName: 'ACP', present: false}, 
        TE:  {domainName: 'TE', present: true},        
      },
      DomainList: 'PKS',
      ModuleCount: 1,
    }
  }

  buildModules = () => {
    let emptyArray = Array.from({ length: this.state.ModuleCount }, ((_, i) =>i+1));
    let moduleArray = emptyArray.map((index) => (
      <ModuleBuilder index={index} key={index} domainList={this.state.PKSDomainList} />
    ));
    return moduleArray;
  }

  addModule = () => {
    let currentPlusOne = this.state.ModuleCount + 1;
    this.setState({ModuleCount: currentPlusOne});
  }

  // select PKS or NRBS
  // delete module button
  // submit button

  render() {
    return (
      <div className='DomainSearch'>
        <h3>Construct Modules</h3>
        <Button onClick={ () => { this.addModule() }}> Add Module + </Button>
        <ModuleBuilder index='0' key='0-loading' domainList={this.state.LoadingPKSList} />
        { this.buildModules() }
        <ModuleBuilder key='terminating' domainList={this.state.TerminatingPKSList} />
        }
      </div>
    )
  }
};

export default DomainSearch;