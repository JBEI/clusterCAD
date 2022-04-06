import React from 'react';
import Button from '../components/Button';
import ModuleBuilder from '../components/ModuleBuilder'

class DomainSearch extends React.Component {

  constructor(props) {
    super(props);
    this.state = {
      PKSDomainList: {
        KS:  {domainName: 'KS', present: false},
        AT:  {domainName: 'AT', present: false},
        DH:  {domainName: 'DH', present: false},
        ER:  {domainName: 'ER', present: false},
        KR:  {domainName: 'KR', present: false},
        ACP: {domainName: 'ACP', present: false}, 
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
      ModuleArray: [],
      DomainList: 'PKS',
      ModuleCount: 0, // does not include loading & terminating
    }
  }

  buildModules = () => {
    let emptyArray = Array.from({ length: this.state.ModuleCount }, ((_, i) =>i+1));
    let moduleArray = emptyArray.map((index) => (
      <ModuleBuilder index={index} key={index} domainList={this.state.PKSDomainList} type="extending" />
    ));
    // this.setState({ModuleArray: moduleArray});
    return moduleArray;
  }

  addModule = () => {
    let currentPlusOne = this.state.ModuleCount + 1;
    this.setState({ModuleCount: currentPlusOne});
  }

  deleteModule = (moduleIndex) => {
    let currentModules = this.state.ModuleArray;
    let splicedArray = currentModules.splice(moduleIndex, 1);
    this.setState({ModuleArray: splicedArray});
    // return currentModules;
  }

  // select PKS or NRBS
  // delete module button
  // submit button

  render() {
    return (
      <div className='DomainSearch'>
        <h3>Construct Modules</h3>
        <Button onClick={ () => { this.addModule() }}> Add Module + </Button>
        <div className="ModuleListWrapper">
          <ModuleBuilder index='0' key='loading' domainList={this.state.LoadingPKSList} type="loading" />
          { this.buildModules() }
          <ModuleBuilder index ={this.state.ModuleCount + 2} key='terminating' domainList={this.state.TerminatingPKSList} type="terminating" />
        </div>
      </div>
    )
  }
};

export default DomainSearch;