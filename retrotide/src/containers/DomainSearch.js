import React from 'react';
import Button from '../components/Button';
import ModuleBuilder from '../components/ModuleBuilder'

class DomainSearch extends React.Component {

  constructor(props) {
    super(props);
    this.state = {
      PKSDomainList: {
        KS:  {domainName: 'KS', present: false},
        AT:  {domainName: 'AT', present: true},
        DH:  {domainName: 'DH', present: false},
        ER:  {domainName: 'ER', present: false},
        KR:  {domainName: 'KR', present: false},
        ACP: {domainName: 'ACP', present: true}, 
      },
      LoadingPKSList: {
        KS:  {domainName: 'KS', present: false},
        AT:  {domainName: 'AT', present: true},
        ACP: {domainName: 'ACP', present: true},       
      },
      TerminatingPKSList: {
        KS:  {domainName: 'KS', present: false},
        AT:  {domainName: 'AT', present: true},
        DH:  {domainName: 'DH', present: false},
        ER:  {domainName: 'ER', present: false},
        KR:  {domainName: 'KR', present: false},
        ACP: {domainName: 'ACP', present: true}, 
        TE:  {domainName: 'TE', present: true},        
      },
      ModuleArray: [],
      DomainList: 'PKS',
    }
  }

  buildModules = (newLength) => {
    let emptyArray = Array.from({ length: newLength }, ((_, i) =>i+1));
    this.setState({ModuleArray: emptyArray});
    console.log("built modules and memoized");
  }

  addModule = () => {
    let currentPlusOne = this.state.ModuleArray.length + 1;
    this.buildModules(currentPlusOne);
  }

  deleteModule = (moduleIndex) => {
    console.log("module to delete has index " + moduleIndex);
    let currentModules = this.state.ModuleArray;
    console.log("modules look like " + currentModules);
    // note that splice() is a *mutator*
    currentModules.splice(moduleIndex-1, 1);
    console.log("modules look like " + currentModules);
    this.setState({
      ModuleArray: currentModules,
    });
  }

  // select PKS or NRBS
  // delete module button
  // submit button
  // key as index is antipattern :(

  render() {
    return (
      <div className='DomainSearch'>
        <h3>Construct Modules</h3>
        <Button onClick={ () => { this.addModule() }}> Add Module + </Button>
        <div className="ModuleListWrapper">
          <ModuleBuilder index='0' key='loading' domainList={this.state.LoadingPKSList} type="loading" />
          { this.state.ModuleArray.map((index) => (
              <ModuleBuilder index={index} key={index} deleteFunction={this.deleteModule} domainList={this.state.PKSDomainList} type="extending" />
            ))
          }
          <ModuleBuilder index ={this.state.ModuleArray.length + 1} key='terminating' domainList={this.state.TerminatingPKSList} type="terminating" />
        </div>
      </div>
    )
  }
};

export default DomainSearch;