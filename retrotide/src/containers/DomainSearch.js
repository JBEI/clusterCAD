import React from 'react';
import { connect } from 'react-redux';
import { beginDomainSearch } from '../redux/actions/actions';
import {clusterCADDomainSearch} from '../redux/middleware/api';
import Button from '../components/Button';
import ModuleBuilder from '../components/ModuleBuilder'

const mapDispatchToProps = dispatch => {
  return {
    beginDomainSearch: modules => dispatch(beginDomainSearch(modules)),
    dispatch,
  }
};

class DomainSearch extends React.Component {

  constructor(props) {
    super(props);

    const LoadingPKSList = {
      KS:  {domainName: 'KS', present: true},
      AT:  {domainName: 'AT', present: true},
      ACP: {domainName: 'ACP', present: true},       
    };

    const TerminatingPKSList = {
      KS:  {domainName: 'KS', present: true},
      AT:  {domainName: 'AT', present: true},
      DH:  {domainName: 'DH', present: false},
      ER:  {domainName: 'ER', present: false},
      KR:  {domainName: 'KR', present: false},
      ACP: {domainName: 'ACP', present: true}, 
      TE:  {domainName: 'TE', present: true},        
    };

    const ButtonList = {
      KS: {domainName: 'KS', domains: ['KS']},
      KR: {domainName: 'KR', domains: ['KR']},
      DH_KR: {domainName: 'DH-KR', domains: ['DH', 'KR']},
      DH_ER_KR: {domainName: 'DH-ER-KR', domains: ['DH', 'ER', 'KR']},
    }

    this.state = {
      ModuleArray: [],
      LoadingModule: {
        key: 'loading-0',
        type: 'loading',
        domainList: LoadingPKSList,
      },
      TerminatingModule: {
        key: 'terminating-n',
        type: 'terminating',
        domainList: TerminatingPKSList,
        buttonList: ButtonList,
      },
      DomainList: 'PKS',
    }
  }

  buildModules = (newLength) => {
    let PKSDomainList = {
      KS:  {domainName: 'KS', present: true},
      AT:  {domainName: 'AT', present: true},
      DH:  {domainName: 'DH', present: false},
      ER:  {domainName: 'ER', present: false},
      KR:  {domainName: 'KR', present: false},
      ACP: {domainName: 'ACP', present: true}, 
    };
    let ButtonList = {
      KS: {domainName: 'KS', domains: ['KS']},
      KR: {domainName: 'KR', domains: ['KR']},
      DH_KR: {domainName: 'DH-KR', domains: ['DH', 'KR']},
      DH_ER_KR: {domainName: 'DH-ER-KR', domains: ['DH', 'ER', 'KR']}
    };
    let defaultArray = Array.from(
      { length: newLength },
      (x, index) => ({
        index: index,
        key: 'extending-' + index,
        domainList: PKSDomainList,
        buttonList: ButtonList,
        type: 'extending',
      })
    );
    this.setState({ ModuleArray: defaultArray });
  }

  addModule = () => {
    let currentPlusOne = this.state.ModuleArray.length + 1;
    this.buildModules(currentPlusOne);
  }

  parseModuleObject = (module, index) => {
    let {key, type, domainList, buttonList} = module;
    return (
      <ModuleBuilder 
        index={index} 
        key={key}
        id={key}
        deleteFunction={this.deleteModule} 
        domainList={domainList}
        buttonList={buttonList}
        type={type} />
    );
  }

  deleteModule = (moduleKey) => {
    let currentModules = this.state.ModuleArray;
    let moduleIndex = currentModules.findIndex((modObj) => modObj.key === moduleKey);
    if(moduleIndex >= 0) {
      currentModules.splice(moduleIndex, 1);
    }
    this.setState({
      ModuleArray: currentModules,
    });
  }

  submitSearch = () => {
    let loading = this.state.LoadingModule;
    let terminating = this.state.TerminatingModule;
    let currentModules = this.state.ModuleArray;
    // push and unshift are mutators
    // this works fine even if ModuleArray has length 0
    currentModules.unshift(loading);
    currentModules.push(terminating);
    // dispatch action
    this.props.beginDomainSearch(currentModules);
    // call async function. is this where this goes?
    clusterCADDomainSearch(currentModules);
  }

  render() {
    const ExtendingArray = this.state.ModuleArray;
    return (
      <div className='DomainSearch'>
        <h3>Construct Modules</h3>
        <Button onClick={() => { this.addModule() }}> Add Module + </Button>
        <Button onClick={() => { this.submitSearch() }} className="submit"> Submit </Button>
        <div className="ModuleListWrapper">
          { this.parseModuleObject(this.state.LoadingModule, -1) }
          { (ExtendingArray && ExtendingArray.length > 0) &&
            ExtendingArray.map((exModule, index) => {
              return (
                this.parseModuleObject(exModule, index)
              )
            }) 
          }
          { this.parseModuleObject(this.state.TerminatingModule, ExtendingArray.length) }
        </div>
      </div>
    )
  }
};

export default connect(null, mapDispatchToProps)(DomainSearch);