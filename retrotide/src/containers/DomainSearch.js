import React from 'react';
import { connect } from 'react-redux';
import { addModule, deleteModule } from '../redux/actions/actions';
import Button from '../components/Button';
import ModuleBuilder from '../components/ModuleBuilder'

const mapDispatchToProps = dispatch => {
  return {
    addModule: moduleList => dispatch(addModule(moduleList)),
    deleteModule: moduleList => dispatch(deleteModule(moduleList)),
    dispatch,
  }
};

class DomainSearch extends React.Component {

  constructor(props) {
    super(props);
    
    const PKSDomainList = {
      KS:  {domainName: 'KS', present: false},
      AT:  {domainName: 'AT', present: true},
      DH:  {domainName: 'DH', present: false},
      ER:  {domainName: 'ER', present: false},
      KR:  {domainName: 'KR', present: false},
      ACP: {domainName: 'ACP', present: true}, 
    };

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
      LoadingModule: {
        index: 0,
        key: 'loading-0',
        type: 'loading',
        domainList: LoadingPKSList,
      },
      // ModuleArray: this.props.ModuleArray,
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
    let emptyArray = Array.from({ length: newLength }, ((_, i) =>i+1));
    emptyArray.map();
    this.setState({ModuleArray: emptyArray});
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
    console.log("modules look like " + currentModules);
    // note that splice() is a *mutator*
    currentModules.splice(moduleIndex-1, 1);
    console.log("modules look like " + currentModules);
    this.setState({
      ModuleArray: currentModules,
    });
  }

  render() {
    return (
      <div className='DomainSearch'>
        <h3>Construct Modules</h3>
        <Button onClick={ () => { this.addModule() }}> Add Module + </Button>
        <div className="ModuleListWrapper">
          { this.parseModuleObject(this.state.LoadingModule) }

          { this.parseModuleObject(this.state.TerminatingModule) }
        </div>
      </div>
    )
  }
};

export default connect(null, mapDispatchToProps) (DomainSearch);