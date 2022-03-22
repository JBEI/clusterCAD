import React from 'react';
import Button from '../components/Button';

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

  getAllDomains = () => {
    let allDomains = [];
    for (var DomainObject in this.state.PKSDomainList) {
      allDomains.push(this.state.PKSDomainList[DomainObject]);
    }
    return allDomains;
  };

  getPresentDomains = () => {
    let presentDomains = [];
    for (var DomainObject in this.state.PKSDomainList) {
      if (this.state.PKSDomainList[DomainObject].present) {
        presentDomains.push(this.state.PKSDomainList[DomainObject]);
      }
    }
    return presentDomains;
  };

  // getPresentDomainsInOrder

  insertDomain = Domain => {
    let selectedDomain = this.state.PKSDomainList[Domain];

    if(selectedDomain.present) {
      console.log("ERR that domain is already selected " + selectedDomain.domainName);
      return -1;
    } else {
      // we'll need some logic here to add multiple nodes depending on selected name
      let insertDomain = {
        domainName: Domain,
        present: true,
      }
      let updatedDomainList = {
        ...this.state.PKSDomainList,
        [Domain]: insertDomain,
      };
      this.setState({PKSDomainList: updatedDomainList});      
    }
  }

  deleteDomain = Domain => {
    let deleteDomain = {
      domainName: Domain,
      present: false,
    }
    let updatedDomainList = {
      ...this.state.PKSDomainList,
      [Domain]: deleteDomain,
    };
    this.setState({PKSDomainList: updatedDomainList});
  }

  render() {
    return (
      <div className='DomainSearch'>
        <h3>Construct Module</h3> 
        <div className="DomainToolbox">
          <div className="DomainButtonList">
            {this.getAllDomains().map((DomainButton, index) => (
              <Button className='addDomainButton' key={index} onClick={ () => {this.insertDomain(DomainButton.domainName)} }>
                {DomainButton.domainName}
              </Button>
              ))
            }
          </div>
          <div className="DomainSandbox">
            {this.getPresentDomains().map((DomainDiv, index) => (
                <div key={DomainDiv.domainName + index}>
                  <div className={"Domain " + DomainDiv.domainName}>
                    {DomainDiv.domainName}
                    <div className="handle"> -> </div>
                    <div className="deleteIcon" onClick={ () => {this.deleteDomain(DomainDiv.domainName)} }> X </div>
                  </div>
                </div>
              ))
            }
          </div>
        </div>
      </div>
    )
  }
};

export default DomainSearch;