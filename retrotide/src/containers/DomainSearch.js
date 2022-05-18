import React from 'react';
import Button from '../components/Button';
import ModuleBuilder from '../components/ModuleBuilder'

class DomainSearch extends React.Component {

  constructor(props) {
    super(props);

    const LoadingPKSList = {
      KS:  {domainName: 'KS', present: false},
      AT:  {domainName: 'AT', present: true},
      ACP: {domainName: 'ACP', present: true},       
    };

    const TerminatingPKSList = {
      KS:  {domainName: 'KS', present: false},
      AT:  {domainName: 'AT', present: true},
      DH:  {domainName: 'DH', present: false},
      ER:  {domainName: 'ER', present: false},
      KR:  {domainName: 'KR', present: false},
      ACP: {domainName: 'ACP', present: true}, 
      TE:  {domainName: 'TE', present: true},        
    };

    this.state = {
      ModuleArray: [],
      LoadingModule: {
        index: 0,
        key: 'loading-0',
        type: 'loading',
        domainList: LoadingPKSList,
      },
      TerminatingModule: {
        index: 1,
        key: 'terminating-n',
        type: 'terminating',
        domainList: TerminatingPKSList,
      },
      DomainList: 'PKS',
    }
  }

  buildModules = (newLength) => {
    let PKSDomainList = {
      KS:  {domainName: 'KS', present: false},
      AT:  {domainName: 'AT', present: true},
      DH:  {domainName: 'DH', present: false},
      ER:  {domainName: 'ER', present: false},
      KR:  {domainName: 'KR', present: false},
      ACP: {domainName: 'ACP', present: true}, 
    };
    let defaultArray = Array.from(
      { length: newLength },
      (x, index) => ({
        index: index,
        key: 'extending-' + index,
        domainList: PKSDomainList,
        type: 'extending',
      })
    );
    console.log(defaultArray);
    this.setState({ ModuleArray: defaultArray });
    console.log("built modules and memoized");
  }

  addModule = () => {
    let currentPlusOne = this.state.ModuleArray.length + 1;
    this.buildModules(currentPlusOne);
  }

  parseModuleObject = (module) => {
    let {index, key, type, domainList} = module;
    return (
      <ModuleBuilder 
        index={index} 
        key={key} 
        deleteFunction={this.deleteModule} 
        domainList={domainList} 
        type={type} />
    );
  }

  deleteModule = (moduleIndex) => {
    console.log("module to delete has index " + moduleIndex);
    let currentModules = this.state.ModuleArray;
    currentModules.splice(moduleIndex, 1);
    this.setState({
      ModuleArray: currentModules,
    });
  }

  render() {
    const ExtendingArray = this.state.ModuleArray;
    return (
      <div className='DomainSearch'>
        <h3>Construct Modules</h3>
        <Button onClick={ () => { this.addModule() }}> Add Module + </Button>
        <div className="ModuleListWrapper">
          { this.parseModuleObject(this.state.LoadingModule) }
          { (ExtendingArray && ExtendingArray.length > 0) &&
            ExtendingArray.map((exModule) => {
              return (
                this.parseModuleObject(exModule)
              )
            }) 
          }
          { this.parseModuleObject(this.state.TerminatingModule) }
        </div>
      </div>
    )
  }
};

export default DomainSearch;